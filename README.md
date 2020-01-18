# webBLAST
Stateless Flask app for running simple BLAST queries.

## What is this?
This is a simple flask app which mostly replicates the functionality of the command line BLAST interface. It's currently a stateless interface, so you can run simple queries and download the results. Eventually, I think I will add user profiles so that you can look at your history and decide which results to keep vs delete.

The functionality is implemented based on the page [here](https://open.oregonstate.education/computationalbiology/chapter/command-line-blast/)

## Why do we have it?

## How does it work?

## How do I operate it?

`gunicorn app:app -b 127.0.0.1:8000 --daemon`

## Where does it live?